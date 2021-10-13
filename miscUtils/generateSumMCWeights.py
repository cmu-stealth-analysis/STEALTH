#!/usr/bin/env python

from __future__ import print_function, division

import subprocess

n_subsamples = {
    "MC_HighHTQCD16": 7,
    "MC_HighHTQCD17": 8,
    "MC_HighHTQCD18": 8,
    "MC_GJetHT16": 5,
    "MC_GJetHT17": 5,
    "MC_GJetHT18": 5
}

input_file_list_and_output_details = {}
for year_last_two_digits in [16, 17, 18]:
    year = 2000 + year_last_two_digits
    input_file_list_and_output_details["MC_DiPhotonJets_{y}".format(y=year)] = ("fileLists/inputFileList_MC_DiPhotonJets_{y}.txt".format(y=year), "xSecLumiInfo", "sumMCWeights_DiPhotonJets_{y}.json".format(y=year))
    for MCBKGDatasetID in ["HighHTQCD", "GJetHT"]:
        for index_subsample in range(1, 1+n_subsamples["MC_{did}{y2}".format(did=MCBKGDatasetID, y2=year_last_two_digits)]):
            input_file_list_and_output_details["MC_{did}{y2}_{i}".format(did=MCBKGDatasetID, y2=year_last_two_digits, i=index_subsample)] = ("fileLists/inputFileList_MC_{did}{i}_{y}.txt".format(did=MCBKGDatasetID, i=index_subsample, y=year), "xSecLumiInfo", "sumMCWeights_{did}_{y}_{i}.json".format(did=MCBKGDatasetID, y=year, i=index_subsample))

for input_file_list_and_output_details_key in input_file_list_and_output_details:
    input_file_list, output_folder, output_file_path = input_file_list_and_output_details[input_file_list_and_output_details_key]
    print("Running getSumMCWeights for input_file_list: {i}, output_folder: {o}, output_file_path: {op}".format(i=input_file_list, o=output_folder, op=output_file_path))
    command_to_call = "./miscUtils/miscScripts/bin/getSumMCWeights inputPathsFiles={i} outputFolder={o} outputFileName={of}".format(i=input_file_list, o=output_folder, of=output_file_path)
    subprocess.check_call(command_to_call, shell=True, executable="/bin/bash")
