#!/usr/bin/env python

from __future__ import print_function, division

import subprocess, os

dataset_names = ["DiPhotonJetsBox", "GJetHT", "HighHTQCD", "data"]
dataset_name_contains_y2 = {
    "GJetHT": True,
    "DiPhotonJetsBox": False,
    "HighHTQCD": True,
    "data": False
}
dataset_is_MC = {
    "GJetHT": True,
    "DiPhotonJetsBox": True,
    "HighHTQCD": True,
    "data": False
}
selections = ["signal", "signal_loose"]
input_folder_with_prefix = "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined_DoublePhoton_MCBkg_lowST/"
minST_nJetsDistributions = 1000.0

input_file_paths_with_prefix_map = {}
for year in [2016, 2017, 2018]:
    year_last_two_digits = year-2000
    input_file_paths_with_prefix_map[str(year)] = {}
    for selection in selections:
        input_file_paths_with_prefix_map[str(year)][selection] = {}
        for dataset_name in dataset_names:
            y2_string = ""
            if (dataset_name_contains_y2[dataset_name]): y2_string = str(year_last_two_digits)
            if (dataset_is_MC[dataset_name]):
                input_file_paths_with_prefix_map[str(year)][selection][dataset_name] = "{i}/merged_selection_MC_{dname}{y2s}_noJetSelection_{y}_{s}.root".format(i=input_folder_with_prefix, dname=dataset_name, y2s=y2_string, y=year, s=selection)
            else:
                input_file_paths_with_prefix_map[str(year)][selection][dataset_name] = "{i}/merged_selection_{dname}{y2s}_noJetSelection_{y}_{s}.root".format(i=input_folder_with_prefix, dname=dataset_name, y2s=y2_string, y=year, s=selection)

input_file_paths_with_prefix_map["all"] = {}
for selection in selections:
    input_file_paths_with_prefix_map["all"][selection] = {}
    for dataset_name in dataset_names:
        input_file_paths_with_prefix_map["all"][selection][dataset_name] = "{f16};{f17};{f18}".format(f16=input_file_paths_with_prefix_map["2016"][selection][dataset_name], f17=input_file_paths_with_prefix_map["2017"][selection][dataset_name], f18=input_file_paths_with_prefix_map["2018"][selection][dataset_name])

output_folder = "/uscms/home/tmudholk/nobackup/analysisAreas/lowSTHistograms_MCBkg_lowST"
if not(os.path.isdir(output_folder)): subprocess.check_call("mkdir -p {o}".format(o=output_folder), shell=True, executable="/bin/bash")

subprocess.check_call("cd /uscms/home/tmudholk/private/stealth/STEALTH/miscUtils/miscScripts && make", shell=True, executable="/bin/bash")

for year_string in ["2016", "2017", "2018", "all"]:
    for selection in selections:
        for dataset_name in dataset_names:
            input_file_paths = input_file_paths_with_prefix_map[year_string][selection][dataset_name]
            output_file_name = "{d}_{y}_{s}.root".format(d=dataset_name, y=year_string, s=selection)
            make_histograms_command = "cd /uscms/home/tmudholk/private/stealth/STEALTH &&"
            make_histograms_command += " ./miscUtils/miscScripts/bin/makeHistogramsLowST"
            make_histograms_command += " \"inputFilePaths={i}\"".format(i=input_file_paths)
            make_histograms_command += " \"outputFolder={o}\"".format(o=output_folder)
            make_histograms_command += " \"outputFileName={o}\"".format(o=output_file_name)
            if (minST_nJetsDistributions > 0.):
                make_histograms_command += " nJetsDistributionsMinST={s:.2f}".format(s=minST_nJetsDistributions)
            if (dataset_is_MC[dataset_name]):
                make_histograms_command += " \"useMCWeights=true\""
            else:
                make_histograms_command += " \"useMCWeights=false\""
            subprocess.check_call(make_histograms_command, shell=True, executable="/bin/bash")

print("All Done!")
