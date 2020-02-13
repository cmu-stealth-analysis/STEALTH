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
    "fileLists/inputFileList_data_SinglePhoton_2016_ntuplizedDec2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016B-17Jul2018_ver2-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016C-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016D-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016E-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016F-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016G-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016H-17Jul2018-v1_singlePhoton_ntuplizedDec2019"
    ],
    "fileLists/inputFileList_data_JetHT_2016_ntuplizedDec2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016B-17Jul2018_ver2-v2_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016C-17Jul2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016D-17Jul2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016E-17Jul2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016F-17Jul2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016G-17Jul2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016H-17Jul2018-v1_JetHT_ntuplizedDec2019"
    ],
    "fileLists/inputFileList_data_DoubleEG_2017_ntuplizedOct2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017B-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017C-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017D-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017E-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017F-31Mar2018-v1_ntuplizedOct2019"
    ],
    "fileLists/inputFileList_data_SinglePhoton_2017_ntuplizedDec2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017B-31Mar2018-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017C-31Mar2018-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017D-31Mar2018-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017E-31Mar2018-v1_singlePhoton_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017F-31Mar2018-v1_singlePhoton_ntuplizedDec2019"
    ],
    "fileLists/inputFileList_data_JetHT_2017_ntuplizedDec2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017B-31Mar2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017C-31Mar2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017D-31Mar2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017E-31Mar2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017F-31Mar2018-v1_JetHT_ntuplizedDec2019",
    ],
    "fileLists/inputFileList_data_EGamma_2018_ntuplizedOct2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018A-17Sep2018-v2_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018B-17Sep2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018C-17Sep2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018D-22Jan2019-v2_ntuplizedOct2019"
    ],
    "fileLists/inputFileList_data_JetHT_2018_ntuplizedDec2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018A-17Sep2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018B-17Sep2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018C-17Sep2018-v1_JetHT_ntuplizedDec2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018D-PromptReco-v2_JetHT_ntuplizedDec2019",
    ],
    "fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_stealth_T5Wg_producedAug19"
    ],
    "fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_stealth_T6Wg_producedAug19"
    ],
    "fileLists/inputFileList_MC_Fall17_MC_DoubleEMEnrichedQCD.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_EMEnrichedQCD_Pt-30to40_producedAug19/",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_EMEnrichedQCD_Pt-30toInf_producedAug19/",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_EMEnrichedQCD_Pt-40toInf_producedAug19/"
    ]
}

for targetFile in eosTargets.keys():
    outputFile = open(targetFile, "w")
    for folder in eosTargets[targetFile]:
        print("Getting list of files in folder: {f}".format(f=folder))
        root_files_list_generator = tmEOSUtils.generate_list_of_files_in_eos_path(eos_path=folder, appendPrefix=True, vetoPattern="failed", restrictToROOTFiles=True, fetchSizeInfo=False)
        for root_file in root_files_list_generator:
            outputFile.write("{rF}\n".format(rF=root_file))
    outputFile.close()
