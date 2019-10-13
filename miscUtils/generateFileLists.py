#!/usr/bin/env python

from __future__ import print_function, division

import tmEOSUtils

eosTargets = {
    "fileLists/inputFileList_data_DoubleEG_2016_ntuplizedOct2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2016B-17Jul2018_ver2-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2016C-17Jul2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2016D-17Jul2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2016E-17Jul2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2016F-17Jul2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2016G-17Jul2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2016H-17Jul2018-v1_ntuplizedOct2019"
    ],
    "fileLists/inputFileList_data_DoubleEG_2017_ntuplizedOct2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017B-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017C-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017D-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017E-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017F-31Mar2018-v1_ntuplizedOct2019"
    ],
    "fileLists/inputFileList_data_EGamma_2018_ntuplizedOct2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018A-17Sep2018-v2_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018B-17Sep2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018C-17Sep2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018D-22Jan2019-v2_ntuplizedOct2019"
    ],
    "fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_stealth_T5Wg_producedAug19"
    ]
}

for targetFile in eosTargets.keys():
    outputFile = open(targetFile, "w")
    for folder in eosTargets[targetFile]:
        print("Getting list of files in folder: {f}".format(f=folder))
        root_files_list_generator = tmEOSUtils.generate_list_of_root_files_in_eos_path(eos_path=folder, appendPrefix=True, vetoPattern="failed")
        for root_file in root_files_list_generator:
            outputFile.write("{rF}\n".format(rF=root_file))
    outputFile.close()
