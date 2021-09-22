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
    # "fileLists/inputFileList_data_SinglePhoton_2016_ntuplizedDec2019.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016B-17Jul2018_ver2-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016C-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016D-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016E-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016F-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016G-17Jul2018-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016H-17Jul2018-v1_singlePhoton_ntuplizedDec2019"
    # ],
    "fileLists/inputFileList_data_SinglePhoton_2016_ntuplizedFeb2021.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2016B-17Jul2018_ver2-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2016C-17Jul2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2016D-17Jul2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2016E-17Jul2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2016F-17Jul2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2016G-17Jul2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2016H-17Jul2018-v1_ntuplizedFeb2021"
    ],
    # "fileLists/inputFileList_data_JetHT_2016_ntuplizedDec2019.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016B-17Jul2018_ver2-v2_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016C-17Jul2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016D-17Jul2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016E-17Jul2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016F-17Jul2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016G-17Jul2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2016H-17Jul2018-v1_JetHT_ntuplizedDec2019"
    # ],
    "fileLists/inputFileList_data_DoubleEG_2017_ntuplizedOct2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017B-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017C-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017D-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017E-31Mar2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/data_Run2017F-31Mar2018-v1_ntuplizedOct2019"
    ],
    # "fileLists/inputFileList_data_SinglePhoton_2017_ntuplizedDec2019.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017B-31Mar2018-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017C-31Mar2018-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017D-31Mar2018-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017E-31Mar2018-v1_singlePhoton_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017F-31Mar2018-v1_singlePhoton_ntuplizedDec2019"
    # ],
    "fileLists/inputFileList_data_SinglePhoton_2017_ntuplizedFeb2021.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2017B-31Mar2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2017C-31Mar2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2017D-31Mar2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2017E-31Mar2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_SinglePhoton_Run2017F-31Mar2018-v1_ntuplizedFeb2021"
    ],
    # "fileLists/inputFileList_data_JetHT_2017_ntuplizedDec2019.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017B-31Mar2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017C-31Mar2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017D-31Mar2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017E-31Mar2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2017F-31Mar2018-v1_JetHT_ntuplizedDec2019",
    # ],
    "fileLists/inputFileList_data_EGamma_2018_ntuplizedOct2019.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018A-17Sep2018-v2_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018B-17Sep2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018C-17Sep2018-v1_ntuplizedOct2019",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018D-22Jan2019-v2_ntuplizedOct2019"
    ],
    "fileLists/inputFileList_data_EGamma_2018_ntuplizedFeb2021.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_EGamma_Run2018A-17Sep2018-v2_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_EGamma_Run2018B-17Sep2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_EGamma_Run2018C-17Sep2018-v1_ntuplizedFeb2021",
        "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_EGamma_Run2018D-22Jan2019-v2_ntuplizedFeb2021"
    ],
    # "fileLists/inputFileList_data_JetHT_2018_ntuplizedDec2019.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018A-17Sep2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018B-17Sep2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018C-17Sep2018-v1_JetHT_ntuplizedDec2019",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/data_Run2018D-PromptReco-v2_JetHT_ntuplizedDec2019",
    # ],
    "fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_stealth_T5Wg_producedAug19"
    ],
    "fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt": [
        "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_stealth_T6Wg_producedAug19"
    ],
    # "fileLists/inputFileList_MC_Fall17_MC_DoubleEMEnrichedQCD.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_EMEnrichedQCD_Pt-30to40_producedAug19/",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_EMEnrichedQCD_Pt-30toInf_producedAug19/",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_EMEnrichedQCD_Pt-40toInf_producedAug19/"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_DoubleEMEnrichedQCD1.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_EMEnrichedQCD_Pt-30to40_producedAug19/"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_DoubleEMEnrichedQCD2.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_EMEnrichedQCD_Pt-30toInf_producedAug19/"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_DoubleEMEnrichedQCD3.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with9413/MC_Fall17_EMEnrichedQCD_Pt-40toInf_producedAug19/"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_QCD.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT300to500",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT500to700",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT700to1000",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT1000to1500",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT1500to2000",
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT2000toInf"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_QCD1.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT300to500"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_QCD2.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT500to700"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_QCD3.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT700to1000"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_QCD4.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT1000to1500"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_QCD5.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT1500to2000"
    # ],
    # "fileLists/inputFileList_MC_Fall17_MC_QCD6.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_QCD_HT2000toInf"
    # ],
    # "fileLists/inputFileList_MC_Summer16_QCD1.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_QCD_HT300to500"
    # ],
    # "fileLists/inputFileList_MC_Summer16_QCD2.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_QCD_HT500to700"
    # ],
    # "fileLists/inputFileList_MC_Summer16_QCD3.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_QCD_HT700to1000"
    # ],
    # "fileLists/inputFileList_MC_Summer16_QCD4.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_QCD_HT1000to1500"
    # ],
    # "fileLists/inputFileList_MC_Summer16_QCD5.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_QCD_HT1500to2000"
    # ],
    # "fileLists/inputFileList_MC_Summer16_QCD6.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_QCD_HT2000toInf"
    # ],
    # "fileLists/inputFileList_MC_Spring18_QCD1.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_QCD_HT300to500"
    # ],
    # "fileLists/inputFileList_MC_Spring18_QCD2.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_QCD_HT500to700"
    # ],
    # "fileLists/inputFileList_MC_Spring18_QCD3.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_QCD_HT700to1000"
    # ],
    # "fileLists/inputFileList_MC_Spring18_QCD4.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_QCD_HT1000to1500"
    # ],
    # "fileLists/inputFileList_MC_Spring18_QCD5.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_QCD_HT1500to2000"
    # ],
    # "fileLists/inputFileList_MC_Spring18_QCD6.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_QCD_HT2000toInf"
    # ],
    # "fileLists/inputFileList_MC_Summer16_GJet1.txt": [
    #     "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_GJetHT_40_100_producedApr21"
    # ],
    # "fileLists/inputFileList_MC_Summer16_GJet2.txt": [
    #     "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_GJetHT_100_200_producedApr21"
    # ],
    # "fileLists/inputFileList_MC_Summer16_GJet3.txt": [
    #     "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_GJetHT_200_400_producedApr21"
    # ],
    # "fileLists/inputFileList_MC_Summer16_GJet4.txt": [
    #     "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_GJetHT_400_600_producedApr21"
    # ],
    # "fileLists/inputFileList_MC_Summer16_GJet5.txt": [
    #     "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_GJetHT_600_inf_producedApr21"
    # ],
    # "fileLists/inputFileList_MC_Spring18_GJet1.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_GJetHT_100_200_producedDec20"
    # ],
    # "fileLists/inputFileList_MC_Spring18_GJet2.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_GJetHT_200_400_producedDec20"
    # ],
    # "fileLists/inputFileList_MC_Spring18_GJet3.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_GJetHT_400_600_producedDec20"
    # ],
    # "fileLists/inputFileList_MC_Spring18_GJet4.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Spring18_GJetHT_600_producedDec20"
    # ],
    # "fileLists/inputFileList_MC_Fall17_GJet1.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_GJetHT_100_200_producedApr21"
    # ],
    # "fileLists/inputFileList_MC_Fall17_GJet2.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_GJetHT_200_400_producedApr21"
    # ],
    # "fileLists/inputFileList_MC_Fall17_GJet3.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_GJetHT_400_600_producedApr21"
    # ],
    # "fileLists/inputFileList_MC_Fall17_GJet4.txt": [
    #     "/store/user/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_GJetHT_600_inf_producedApr21"
    # ],
    "fileLists/inputFileList_MC_Summer16_hgg.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_Summer16_hgg_producedMay21"
    ],
    "fileLists/inputFileList_MC_Fall17_hgg.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_Fall17_hgg_producedJun20"
    ],
    "fileLists/inputFileList_MC_Autumn18_hgg.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_Autumn18_hgg_producedMay21"
    ],
    "fileLists/inputFileList_MC_EMEnrichedGJet1_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_EMEnrichedGJetPt_2016_1_ntuplizedSep2021",
    ],
    "fileLists/inputFileList_MC_EMEnrichedGJet2_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_EMEnrichedGJetPt_2016_2_ntuplizedSep2021",
    ],
    "fileLists/inputFileList_MC_EMEnrichedGJet3_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_EMEnrichedGJetPt_2016_3_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_EMEnrichedGJet1_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_EMEnrichedGJetPt_2017_1_ntuplizedSep2021",
    ],
    "fileLists/inputFileList_MC_EMEnrichedGJet2_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_EMEnrichedGJetPt_2017_2_ntuplizedSep2021",
    ],
    "fileLists/inputFileList_MC_EMEnrichedGJet3_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_EMEnrichedGJetPt_2017_3_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_EMEnrichedGJet1_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_EMEnrichedGJetPt_2018_1_ntuplizedSep2021",
    ],
    "fileLists/inputFileList_MC_EMEnrichedGJet2_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_EMEnrichedGJetPt_2018_2_ntuplizedSep2021",
    ],
    "fileLists/inputFileList_MC_EMEnrichedGJet3_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_EMEnrichedGJetPt_2018_3_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD1_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2016_1_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD2_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2016_2_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD3_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2016_3_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD4_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2016_4_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD5_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2016_5_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD6_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2016_6_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD7_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2016_7_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD1_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2017_1_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD2_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2017_2_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD3_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2017_3_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD4_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2017_4_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD5_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2017_5_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD6_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2017_6_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD7_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2017_7_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD8_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2017_8_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD1_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2018_1_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD2_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2018_2_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD3_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2018_3_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD4_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2018_4_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD5_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2018_5_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD6_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2018_6_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD7_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2018_7_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_HighHTQCD8_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_HighHTQCD_2018_8_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_DiPhotonJets_2016.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_DiPhotonJets_2016_1_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_DiPhotonJets_2017.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_DiPhotonJets_2017_1_ntuplizedSep2021"
    ],
    "fileLists/inputFileList_MC_DiPhotonJets_2018.txt": [
        "/store/group/lpcsusystealth/stealth2018Ntuples_with10210/MC_DiPhotonJets_2018_1_ntuplizedSep2021"
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
